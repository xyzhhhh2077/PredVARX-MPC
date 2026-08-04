"""Verify the encoded copyBI summary GIF and export a keyframe contact sheet."""

import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageChops, ImageDraw

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results_smpc"
GIF_PATH = RESULTS / "copyBI_smpc_summary_200frames.gif"
SHEET_PATH = RESULTS / "copyBI_smpc_summary_keyframes.png"
REPORT_PATH = RESULTS / "copyBI_smpc_summary_gif_audit.json"
EXPECTED_FRAMES = 200


def main():
    image = Image.open(GIF_PATH)
    count = int(image.n_frames)
    if count != EXPECTED_FRAMES:
        raise AssertionError(f"Encoded frame count {count} != {EXPECTED_FRAMES}")

    indices = [0, count // 3, 2 * count // 3, count - 1]
    frames = []
    durations = []
    for index in range(count):
        image.seek(index)
        durations.append(int(image.info.get("duration", 0)))
        if index in indices:
            frames.append(image.convert("RGB").copy())

    width, height = frames[0].size
    sheet = Image.new("RGB", (2 * width, 2 * height), "white")
    draw = ImageDraw.Draw(sheet)
    for cell, (index, frame) in enumerate(zip(indices, frames)):
        x = (cell % 2) * width
        y = (cell // 2) * height
        sheet.paste(frame, (x, y))
        draw.rectangle((x + 8, y + 8, x + 116, y + 34), fill="white", outline="#6B7280")
        draw.text((x + 15, y + 14), f"frame {index + 1}/{count}", fill="#111827")
    sheet.save(SHEET_PATH)

    arrays = [np.asarray(frame, dtype=np.int16) for frame in frames]
    mad_from_first = [float(np.abs(array - arrays[0]).mean()) for array in arrays]
    differences = []
    for previous, current in zip(frames[:-1], frames[1:]):
        differences.append(ImageChops.difference(previous, current).getbbox() is not None)
    if not all(differences):
        raise AssertionError("At least one sampled keyframe pair is pixel-identical")
    if mad_from_first[-1] < 1.0:
        raise AssertionError("Start-to-end motion energy is too low")

    final = arrays[-1]
    luminance = final.mean(axis=2)
    report = {
        "gif": str(GIF_PATH),
        "frames": count,
        "size_pixels": [width, height],
        "duration_ms_unique": sorted(set(durations)),
        "total_duration_seconds": sum(durations) / 1000.0,
        "loop": int(image.info.get("loop", -1)),
        "file_bytes": GIF_PATH.stat().st_size,
        "keyframe_indices_zero_based": indices,
        "keyframe_mad_from_first": mad_from_first,
        "final_frame_mean_luminance": float(luminance.mean()),
        "final_frame_edge_energy": float(np.abs(np.diff(final.astype(float), axis=1)).mean()),
        "contact_sheet": str(SHEET_PATH),
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
