const pptxgen = require("pptxgenjs");

let pres = new pptxgen();
pres.layout = "LAYOUT_16x9";
pres.author = "XIE Yuezheng";
pres.title = "PredVARX-SMPC 算法演进与修正";

// === Color palette ===
const C = {
  navy: "065A82",
  teal: "1C7293",
  midnight: "21295C",
  white: "FFFFFF",
  offWhite: "F0F4F8",
  lightGray: "94A3B8",
  slate: "64748B",
  red: "DC2626",
  green: "16A34A",
  amber: "D97706",
  bg: "0F172A",
  cardBg: "1E293B",
};

const mkSh = () => ({
  type: "outer", color: "000000", blur: 6, offset: 2, angle: 135, opacity: 0.15,
});

// Helper: dark card with left accent
function accentCard(slide, x, y, w, h, accentColor, title, lines) {
  slide.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: C.white }, shadow: mkSh() });
  slide.addShape(pres.shapes.RECTANGLE, { x, y, w: 0.08, h, fill: { color: accentColor } });
  let items = [
    { text: title, options: { fontSize: 14, bold: true, color: accentColor, breakLine: true } },
    { text: " ", options: { fontSize: 4, breakLine: true } },
  ];
  for (let i = 0; i < lines.length; i++) {
    items.push({
      text: lines[i],
      options: { fontSize: 12, color: C.midnight, breakLine: i < lines.length - 1 },
    });
  }
  slide.addText(items, { x: x + 0.25, y: y + 0.15, w: w - 0.45, h: h - 0.25 });
}

// Helper: number circle
function numCircle(slide, x, y, num, color) {
  slide.addShape(pres.shapes.OVAL, { x, y, w: 0.55, h: 0.55, fill: { color } });
  slide.addText(String(num), {
    x, y, w: 0.55, h: 0.55,
    fontSize: 20, bold: true, color: C.white, align: "center", valign: "middle",
  });
}

// Helper: slide header
function slideHeader(slide, title, accentColor) {
  slide.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: 10, h: 0.06, fill: { color: accentColor } });
  slide.addText(title, {
    x: 0.6, y: 0.2, w: 8.8, h: 0.55,
    fontSize: 28, bold: true, color: C.midnight,
  });
}

// ============================================================
// SLIDE 1: Title
// ============================================================
let s1 = pres.addSlide();
s1.background = { color: C.bg };
s1.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: 10, h: 0.06, fill: { color: C.teal } });

s1.addText("PredVARX-SMPC 算法演进与修正", {
  x: 0.8, y: 1.0, w: 8.4, h: 0.9,
  fontSize: 38, fontFace: "Arial", bold: true, color: C.white, align: "left",
});
s1.addText("从离线大数据累积辨识到在线小样本滑动窗口 SMPC", {
  x: 0.8, y: 1.9, w: 8.4, h: 0.6,
  fontSize: 20, fontFace: "Arial", color: C.teal, align: "left",
});

// Four key numbers
let nums = [
  { label: "离线数据", old: "5000步", now: "200步", pct: "−96%" },
  { label: "MAE", old: "0.469", now: "0.395", pct: "+16%" },
  { label: "越界", old: "1", now: "0", pct: "−100%" },
  { label: "SMPC绑定", old: "0%", now: "2.4%", pct: "✓" },
];

for (let i = 0; i < 4; i++) {
  let cx = 0.8 + i * 2.2;
  s1.addShape(pres.shapes.RECTANGLE, {
    x: cx, y: 3.0, w: 2.0, h: 1.6,
    fill: { color: C.cardBg }, shadow: mkSh(),
  });
  s1.addText(nums[i].label, {
    x: cx, y: 3.05, w: 2.0, h: 0.35,
    fontSize: 11, color: C.lightGray, align: "center",
  });
  s1.addText(nums[i].pct, {
    x: cx, y: 3.35, w: 2.0, h: 0.55,
    fontSize: 26, bold: true, color: C.green, align: "center",
  });
  s1.addText(nums[i].old + " → " + nums[i].now, {
    x: cx, y: 3.9, w: 2.0, h: 0.35,
    fontSize: 10, color: C.slate, align: "center",
  });
}

s1.addText("谢跃争 · PredVARX-MPC · 2026.07", {
  x: 0.8, y: 5.0, w: 8.4, h: 0.4,
  fontSize: 12, color: C.lightGray,
});

// ============================================================
// SLIDE 2: Algorithm Overview — Full Pipeline
// ============================================================
let s2 = pres.addSlide();
s2.background = { color: C.offWhite };
slideHeader(s2, "算法总体架构：离线辨识 + 在线 SMPC", C.navy);

// Pipeline flow: 4 boxes connected by arrows
let pipeline = [
  { title: "离线数据采集", sub: "强激励开环", color: C.navy },
  { title: "系统辨识", sub: "IVR + VARX", color: C.teal },
  { title: "在线 MPC", sub: "QP 求解", color: C.navy },
  { title: "SMPC 约束", sub: "机会约束收紧", color: C.green },
];

for (let i = 0; i < 4; i++) {
  let cx = 0.5 + i * 2.45;
  s2.addShape(pres.shapes.RECTANGLE, {
    x: cx, y: 1.0, w: 2.1, h: 1.0,
    fill: { color: pipeline[i].color }, shadow: mkSh(),
  });
  s2.addText([
    { text: pipeline[i].title, options: { fontSize: 14, bold: true, color: C.white, breakLine: true } },
    { text: pipeline[i].sub, options: { fontSize: 11, color: C.lightGray } },
  ], { x: cx, y: 1.0, w: 2.1, h: 1.0, align: "center", valign: "middle" });
  if (i < 3) {
    s2.addShape(pres.shapes.LINE, {
      x: cx + 2.15, y: 1.5, w: 0.25, h: 0,
      line: { color: C.lightGray, width: 2 },
    });
  }
}

// Below: two columns explaining offline vs online
// Left: offline
s2.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 2.3, w: 4.3, h: 3.0,
  fill: { color: C.white }, shadow: mkSh(),
});
s2.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 2.3, w: 4.3, h: 0.06, fill: { color: C.navy },
});
s2.addText([
  { text: "离线阶段（Algorithm 1）", options: { fontSize: 15, bold: true, color: C.navy, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "① OLS 剥离输入响应 u → 残差 y_r", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "② 预白化 + IVR 迭代精炼 → 子空间 P̂, R̂", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "③ 从原始 y 提取降维状态 x̂ = R̂'·y_c", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "④ VARX 回归 → Â, B̂（输入同步中心化）", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "⑤ 噪声解耦：Σ_ε（DLV）+ Σ_ē（静态）", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "输出：Â, B̂, P̂, R̂, Σ_ε, Σ_ē, P̄", options: { fontSize: 12, bold: true, color: C.teal } },
], { x: 0.8, y: 2.45, w: 3.8, h: 2.7 });

// Right: online
s2.addShape(pres.shapes.RECTANGLE, {
  x: 5.2, y: 2.3, w: 4.3, h: 3.0,
  fill: { color: C.white }, shadow: mkSh(),
});
s2.addShape(pres.shapes.RECTANGLE, {
  x: 5.2, y: 2.3, w: 4.3, h: 0.06, fill: { color: C.green },
});
s2.addText([
  { text: "在线阶段（每步 k）", options: { fontSize: 15, bold: true, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "① 观测 y_k → 降维状态 x̂_k = R̂'·y_k", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "② 偏差修正：b_k = EMA(y_k − r_k)", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "③ 新息 ε_k → 滑动窗口估计 σ_ε", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "④ 构建 QP（差分预测 + 输入约束）", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "⑤ SMPC 机会约束：μ + z_α·σ ≤ y_max", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "⑥ 每 Nr 步滑动窗口重辨识 + P̄ 更新", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "输出：u_k（控制输入）", options: { fontSize: 12, bold: true, color: C.teal } },
], { x: 5.5, y: 2.45, w: 3.8, h: 2.7 });

// ============================================================
// SLIDE 3: Problem 1 — Too much offline data
// ============================================================
let s3 = pres.addSlide();
s3.background = { color: C.offWhite };
slideHeader(s3, "问题一：离线数据需求过大", C.red);

numCircle(s3, 0.6, 0.9, 1, C.red);

// Before/after comparison
s3.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.6, w: 4.2, h: 2.0,
  fill: { color: C.white }, shadow: mkSh(),
});
s3.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.6, w: 0.08, h: 2.0, fill: { color: C.red },
});
s3.addText([
  { text: "之前", options: { fontSize: 16, bold: true, color: C.red, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "T_off = 5000 ~ 20000 步", options: { fontSize: 22, bold: true, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "需要长时间开环强激励采集", options: { fontSize: 13, color: C.slate, breakLine: true } },
  { text: "实际工程中成本高、条件苛刻", options: { fontSize: 13, color: C.slate } },
], { x: 0.9, y: 1.65, w: 3.6, h: 1.9 });

s3.addShape(pres.shapes.RECTANGLE, {
  x: 5.3, y: 1.6, w: 4.2, h: 2.0,
  fill: { color: C.white }, shadow: mkSh(),
});
s3.addShape(pres.shapes.RECTANGLE, {
  x: 5.3, y: 1.6, w: 0.08, h: 2.0, fill: { color: C.green },
});
s3.addText([
  { text: "修正后", options: { fontSize: 16, bold: true, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "T_off = 200 步", options: { fontSize: 22, bold: true, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "数据需求降低 96%", options: { fontSize: 13, color: C.slate, breakLine: true } },
  { text: "在线滑动窗口持续补充当前工况", options: { fontSize: 13, color: C.slate } },
], { x: 5.7, y: 1.65, w: 3.6, h: 1.9 });

// Key insight
s3.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 3.9, w: 9.0, h: 1.3,
  fill: { color: C.cardBg }, shadow: mkSh(),
});
s3.addText([
  { text: "为什么 200 步就够了？", options: { fontSize: 15, bold: true, color: C.teal, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "• 离线数据只用于初始辨识（IVR 收敛 + VARX 回归），不需要覆盖全部工况", options: { fontSize: 12, color: C.white, breakLine: true } },
  { text: "• 在线阶段每 Nr=30 步触发滑动窗口重辨识，持续用当前数据修正模型", options: { fontSize: 12, color: C.white, breakLine: true } },
  { text: "• 离线 200 步 + 在线滑动窗口 = 混合窗口策略，兼顾初始精度与实时适应", options: { fontSize: 12, color: C.white } },
], { x: 0.8, y: 3.95, w: 8.4, h: 1.2 });

// ============================================================
// SLIDE 4: Problem 2 — Cumulative identification
// ============================================================
let s4 = pres.addSlide();
s4.background = { color: C.offWhite };
slideHeader(s4, "问题二：累积辨识导致模型退化", C.red);

numCircle(s4, 0.6, 0.9, 2, C.red);

// Diagram: cumulative vs sliding
s4.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.6, w: 4.2, h: 2.5,
  fill: { color: C.white }, shadow: mkSh(),
});
s4.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.6, w: 0.08, h: 2.5, fill: { color: C.red },
});
s4.addText([
  { text: "累积辨识（之前）", options: { fontSize: 15, bold: true, color: C.red, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "数据只增不减：", options: { fontSize: 13, bold: true, color: C.midnight, breakLine: true } },
  { text: "[离线 | t=1 | t=2 | ... | t=k]", options: { fontSize: 11, color: C.slate, breakLine: true } },
  { text: "← 全部用于重辨识 →", options: { fontSize: 11, color: C.red, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "问题：", options: { fontSize: 13, bold: true, color: C.midnight, breakLine: true } },
  { text: "• 旧数据包含过时动态信息", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "• 噪声工况变化后旧数据拖累新模型", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "• 数据量持续增长 → 计算量增加", options: { fontSize: 12, color: C.midnight } },
], { x: 0.9, y: 1.65, w: 3.6, h: 2.4 });

s4.addShape(pres.shapes.RECTANGLE, {
  x: 5.3, y: 1.6, w: 4.2, h: 2.5,
  fill: { color: C.white }, shadow: mkSh(),
});
s4.addShape(pres.shapes.RECTANGLE, {
  x: 5.3, y: 1.6, w: 0.08, h: 2.5, fill: { color: C.green },
});
s4.addText([
  { text: "滑动窗口重辨识（修正后）", options: { fontSize: 15, bold: true, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "混合窗口策略：", options: { fontSize: 13, bold: true, color: C.midnight, breakLine: true } },
  { text: "[离线200步 | 最近100步]", options: { fontSize: 11, color: C.green, breakLine: true } },
  { text: "← 固定基础 + 时变窗口 →", options: { fontSize: 11, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "优势：", options: { fontSize: 13, bold: true, color: C.midnight, breakLine: true } },
  { text: "• 离线数据提供稳定基础动态", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "• 在线窗口捕捉当前工况变化", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "• 计算量恒定，不随时间增长", options: { fontSize: 12, color: C.midnight } },
], { x: 5.7, y: 1.65, w: 3.6, h: 2.4 });

// Pbar sync
s4.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 4.4, w: 9.0, h: 0.9,
  fill: { color: C.cardBg }, shadow: mkSh(),
});
s4.addText([
  { text: "同步更新：", options: { fontSize: 14, bold: true, color: C.teal } },
  { text: "每次重辨识后，P̄ = null(P̂') 确保静态噪声基底与新子空间正交，避免信号泄漏到噪声估计中", options: { fontSize: 12, color: C.white } },
], { x: 0.8, y: 4.45, w: 8.4, h: 0.8, valign: "middle" });

// ============================================================
// SLIDE 5: Problem 3 — Δ-MPC SMPC logic
// ============================================================
let s5 = pres.addSlide();
s5.background = { color: C.offWhite };
slideHeader(s5, "问题三：增量 MPC 与 SMPC 的逻辑矛盾", C.red);

numCircle(s5, 0.6, 0.9, 3, C.red);

// Explanation
s5.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.5, w: 9.0, h: 1.8,
  fill: { color: C.white }, shadow: mkSh(),
});
s5.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.5, w: 0.08, h: 1.8, fill: { color: C.amber },
});
s5.addText([
  { text: "增量 MPC（Δ-MPC）的预测形式", options: { fontSize: 15, bold: true, color: C.amber, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "j=1:  Δr₁ = r̃₁ − ŷ_k  （当前输出到目标的增量）", options: { fontSize: 13, color: C.midnight, breakLine: true } },
  { text: "j≥2: Δr_j = r̃_j − r̃_{j-1} （相邻目标的增量）", options: { fontSize: 13, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "QP 只看到\"目标变化了多少\"，看不到\"系统偏离安全边界多少\"", options: { fontSize: 13, bold: true, color: C.red } },
], { x: 0.9, y: 1.55, w: 8.4, h: 1.7 });

// The contradiction
s5.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 3.5, w: 4.2, h: 1.7,
  fill: { color: C.white }, shadow: mkSh(),
});
s5.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 3.5, w: 4.2, h: 0.06, fill: { color: C.red },
});
s5.addText([
  { text: "SMPC 需要绝对预测", options: { fontSize: 14, bold: true, color: C.red, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "机会约束：μ_y(j) + z_α·σ(j) ≤ y_max", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "μ_y 是绝对预测值，需要知道系统状态", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "在输出空间中的绝对位置", options: { fontSize: 12, color: C.midnight } },
], { x: 0.8, y: 3.55, w: 3.7, h: 1.6 });

s5.addShape(pres.shapes.RECTANGLE, {
  x: 5.3, y: 3.5, w: 4.2, h: 1.7,
  fill: { color: C.white }, shadow: mkSh(),
});
s5.addShape(pres.shapes.RECTANGLE, {
  x: 5.3, y: 3.5, w: 4.2, h: 0.06, fill: { color: C.red },
});
s5.addText([
  { text: "Δ-MPC 只有相对增量", options: { fontSize: 14, bold: true, color: C.red, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "差分形式消除了绝对位置信息：", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "Δr 只反映\"目标变了多少\"", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "无法判断\"系统离边界多远\"", options: { fontSize: 12, color: C.midnight } },
], { x: 5.6, y: 3.55, w: 3.7, h: 1.6 });

// ============================================================
// SLIDE 6: Problem 4 — Input mean centering
// ============================================================
let s6 = pres.addSlide();
s6.background = { color: C.offWhite };
slideHeader(s6, "问题四：输入未中心化导致 B̂ 辨识偏差", C.red);

numCircle(s6, 0.6, 0.9, 4, C.red);

// Math explanation
s6.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.5, w: 9.0, h: 1.6,
  fill: { color: C.white }, shadow: mkSh(),
});
s6.addText([
  { text: "输入非零均值的后果", options: { fontSize: 15, bold: true, color: C.amber, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "真实系统：x_{k+1} = A·x_k + B·u_k + w_k", options: { fontSize: 13, color: C.navy, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  { text: "当 mean(u) ≠ 0 时，VARX 回归将 u 的直流分量混入状态转移：", options: { fontSize: 12, color: C.midnight, breakLine: true } },
  { text: "Â 和 B̂ 的估计产生系统性偏差，特别是 B̂ 的方向被扭曲", options: { fontSize: 12, color: C.midnight } },
], { x: 0.8, y: 1.55, w: 8.4, h: 1.5 });

// Fix
s6.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 3.3, w: 9.0, h: 1.8,
  fill: { color: C.cardBg }, shadow: mkSh(),
});
s6.addText([
  { text: "修正：同步零均值中心化", options: { fontSize: 15, bold: true, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "y_c = y_data − mean(y_data, 2)", options: { fontSize: 14, color: C.white, breakLine: true } },
  { text: "u_c = u_data − mean(u_data, 2)   ← 关键：输入也要中心化", options: { fontSize: 14, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 6, breakLine: true } },
  { text: "VARX 回归使用中心化后的 (y_c, u_c)，消除直流偏置对 Â, B̂ 的影响", options: { fontSize: 12, color: C.lightGray, breakLine: true } },
  { text: "在线阶段 u_k 同样减去均值后再用于状态预测", options: { fontSize: 12, color: C.lightGray } },
], { x: 0.8, y: 3.35, w: 8.4, h: 1.7 });

// ============================================================
// SLIDE 7: Current Algorithm — Complete Flow
// ============================================================
let s7 = pres.addSlide();
s7.background = { color: C.offWhite };
slideHeader(s7, "修正后的完整算法流程", C.green);

// Left column: offline
s7.addShape(pres.shapes.RECTANGLE, {
  x: 0.3, y: 1.0, w: 4.6, h: 4.3,
  fill: { color: C.white }, shadow: mkSh(),
});
s7.addShape(pres.shapes.RECTANGLE, {
  x: 0.3, y: 1.0, w: 4.6, h: 0.06, fill: { color: C.navy },
});

let offlineSteps = [
  "① 采集 T_off=200 步开环强激励数据",
  "② OLS 剥离：y_r = y − Φ_u·Ĉ_OLS",
  "③ 预白化：Y* = D^{−1/2}U'·y_r",
  "④ IVR 迭代（10~30步收敛）",
  "     X̂ = Y*_cur'·Ĉ → Â_var → Π 更新",
  "⑤ 反归一化 + QR：P̂, R̂",
  "⑥ 从原始 y 提取：x̂ = R̂'·y_c",
  "⑦ VARX：[A,B]' = (Φ·Φ'+εI)^{−1}·Φ·X_n'",
  "     ★ 输入同步中心化 u_c = u − mean(u)",
  "⑧ 噪声：Σ_ε = cov(x̂−Âx̂−B̂u), Σ_ē = P̄'·resid",
];
s7.addText([
  { text: "离线阶段", options: { fontSize: 15, bold: true, color: C.navy, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  ...offlineSteps.map((s, i) => ({
    text: s,
    options: { fontSize: 10.5, color: s.startsWith("     ") ? C.teal : C.midnight, breakLine: i < offlineSteps.length - 1 },
  })),
], { x: 0.5, y: 1.1, w: 4.2, h: 4.1 });

// Right column: online
s7.addShape(pres.shapes.RECTANGLE, {
  x: 5.1, y: 1.0, w: 4.6, h: 4.3,
  fill: { color: C.white }, shadow: mkSh(),
});
s7.addShape(pres.shapes.RECTANGLE, {
  x: 5.1, y: 1.0, w: 4.6, h: 0.06, fill: { color: C.green },
});

let onlineSteps = [
  "每步 k = 1, 2, ..., T_cl:",
  " ",
  "① 观测：x̂_k = R̂'·y_k",
  "② 偏差修正：b_k = EMA(y−r)",
  "③ 新息：ε_k = x̂_k − (Â·x̂_{k-1}+B̂·u)",
  "     滑动窗口估计 σ_ε",
  "④ 构建 QP（差分预测 + 约束）",
  "⑤ SMPC：μ_i(j)+z_α·√Σ_ii(j) ≤ y_max",
  "⑥ 求解 u_k = quadprog(H,f,A_cc,b_cc)",
  " ",
  "每 Nr=30 步触发重辨识：",
  "  数据 = [离线200 | 最近100步]",
  "  更新 Â,P̂,M_j,N_j,P̄",
];
s7.addText([
  { text: "在线阶段", options: { fontSize: 15, bold: true, color: C.green, breakLine: true } },
  { text: " ", options: { fontSize: 4, breakLine: true } },
  ...onlineSteps.map((s, i) => ({
    text: s,
    options: {
      fontSize: 10.5,
      color: s.startsWith("  ") ? C.teal : (s === " " ? C.white : C.midnight),
      breakLine: i < onlineSteps.length - 1,
    },
  })),
], { x: 5.3, y: 1.1, w: 4.2, h: 4.1 });

// ============================================================
// SLIDE 8: Current Algorithm — Theory
// ============================================================
let s8 = pres.addSlide();
s8.background = { color: C.offWhite };
slideHeader(s8, "修正后的理论框架", C.navy);

// DLV model
s8.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 1.0, w: 9.0, h: 1.3,
  fill: { color: C.cardBg }, shadow: mkSh(),
});
s8.addText([
  { text: "动态潜变量（DLV）模型", options: { fontSize: 14, bold: true, color: C.teal, breakLine: true } },
  { text: "x_{k+1} = A·x_k + B·u_k + w_k    (ℓ维降维状态)    |    y_k = P·x_k + P̄·ē_k    (p维观测)", options: { fontSize: 13, color: C.white } },
], { x: 0.8, y: 1.05, w: 8.4, h: 1.2, valign: "middle" });

// Three theory blocks
let theoryBlocks = [
  {
    title: "子空间辨识（IVR）",
    color: C.navy,
    lines: [
      "工具变量迭代精炼：",
      "M_IVR = W·(X_pred'X_pred)^{-1}·W'",
      "收敛判据：|tr(Σ_x)变化| < tol",
      "P̂'P̂ = I_ℓ（列半正交）",
      "R̂'P̂ = I_ℓ（斜投影对齐）",
    ],
  },
  {
    title: "MPC 代价函数",
    color: C.teal,
    lines: [
      "J = Σ_j ‖μ_y(j)−r̃_j‖²_Q + ‖u‖²_R",
      "差分预测避免输入-噪声相关",
      "H = Σ N_d'QN_d + I⊗R_w",
      "完整代价 = QP部分 + J_const",
      "Warm-start：上一步解平移",
    ],
  },
  {
    title: "SMPC 机会约束",
    color: C.green,
    lines: [
      "协方差传播：",
      "Σ_z(j) = A·Σ_z(j−1)·A' + G·Σ_ε·G'",
      "Σ_y(j) = P·Σ_z(j)·P' + σ_e²·I",
      "约束：μ_i(j) + z_α·√Σ_ii(j) ≤ y_max",
      "两阶段QP：无约束→约束修正",
    ],
  },
];

for (let i = 0; i < 3; i++) {
  let cx = 0.5 + i * 3.15;
  s8.addShape(pres.shapes.RECTANGLE, {
    x: cx, y: 2.5, w: 2.9, h: 2.8,
    fill: { color: C.white }, shadow: mkSh(),
  });
  s8.addShape(pres.shapes.RECTANGLE, {
    x: cx, y: 2.5, w: 2.9, h: 0.06, fill: { color: theoryBlocks[i].color },
  });
  s8.addText([
    { text: theoryBlocks[i].title, options: { fontSize: 14, bold: true, color: theoryBlocks[i].color, breakLine: true } },
    { text: " ", options: { fontSize: 4, breakLine: true } },
    ...theoryBlocks[i].lines.map((l, j) => ({
      text: l,
      options: { fontSize: 10.5, color: C.midnight, breakLine: j < theoryBlocks[i].lines.length - 1 },
    })),
  ], { x: cx + 0.15, y: 2.6, w: 2.6, h: 2.6 });
}

// ============================================================
// SLIDE 9: Results + Summary
// ============================================================
let s9 = pres.addSlide();
s9.background = { color: C.bg };
s9.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: 10, h: 0.06, fill: { color: C.teal } });

s9.addText("实验结果", {
  x: 0.6, y: 0.3, w: 8, h: 0.6,
  fontSize: 32, bold: true, color: C.white,
});

// Result table
s9.addTable([
  [
    { text: "指标", options: { bold: true, color: C.white, fill: { color: C.midnight }, fontSize: 13 } },
    { text: "副本 C（旧方案）", options: { bold: true, color: C.white, fill: { color: C.red }, fontSize: 13 } },
    { text: "副本 S（修正后）", options: { bold: true, color: C.white, fill: { color: C.green }, fontSize: 13 } },
  ],
  [
    { text: "离线数据", options: { fontSize: 12, color: C.lightGray } },
    { text: "5000步 累积", options: { fontSize: 12, color: C.red } },
    { text: "200步 + 滑动窗口", options: { fontSize: 12, color: C.green } },
  ],
  [
    { text: "辨识策略", options: { fontSize: 12, color: C.lightGray } },
    { text: "累积（只增不减）", options: { fontSize: 12, color: C.red } },
    { text: "混合窗口（离线+在线）", options: { fontSize: 12, color: C.green } },
  ],
  [
    { text: "输入中心化", options: { fontSize: 12, color: C.lightGray } },
    { text: "❌ 未中心化", options: { fontSize: 12, color: C.red } },
    { text: "✅ u_c = u − mean(u)", options: { fontSize: 12, color: C.green } },
  ],
  [
    { text: "SMPC 绑定率", options: { fontSize: 12, color: C.lightGray } },
    { text: "0%（y_max=10 未触发）", options: { fontSize: 12, color: C.red } },
    { text: "2.4%（y_max=3.10 真正触发）", options: { fontSize: 12, color: C.green } },
  ],
  [
    { text: "MAE (y₂)", options: { fontSize: 12, color: C.lightGray } },
    { text: "0.469", options: { fontSize: 12, color: C.red } },
    { text: "0.395（+16%）", options: { fontSize: 12, color: C.green } },
  ],
  [
    { text: "越界次数", options: { fontSize: 12, color: C.lightGray } },
    { text: "1", options: { fontSize: 12, color: C.red } },
    { text: "0（−100%）", options: { fontSize: 12, color: C.green } },
  ],
], {
  x: 0.5, y: 1.1, w: 9, h: 2.8,
  border: { pt: 0.5, color: "334155" },
  colW: [2.2, 3.4, 3.4],
});

// Summary message
s9.addShape(pres.shapes.RECTANGLE, {
  x: 0.5, y: 4.2, w: 9.0, h: 1.1,
  fill: { color: C.cardBg }, shadow: mkSh(),
});
s9.addText([
  { text: "核心改进：", options: { fontSize: 14, bold: true, color: C.teal, breakLine: true } },
  { text: "小样本离线 + 滑动窗口在线重辨识 + 输入中心化 + SMPC 真正触发 → 数据效率提升 25×，控制精度 +16%", options: { fontSize: 13, color: C.white } },
], { x: 0.8, y: 4.25, w: 8.4, h: 1.0, valign: "middle" });

// Save
const outPath = "E:/academic_files/phd-learning/PredVAR+MPC/论文笔记/SMPC算法演进.pptx";
pres.writeFile({ fileName: outPath }).then(() => {
  console.log("Saved: " + outPath);
}).catch(err => {
  console.error("Error:", err);
});
