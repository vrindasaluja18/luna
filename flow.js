/* Two animations for the Luna page.

   c-points  — one cloud read three ways: points -> density -> clusters.
   c-merge   — two peaks merging as a threshold descends past the valley.

   Plain canvas 2D. No libraries, no network. Both respect
   prefers-reduced-motion by drawing a single informative frame. */

(function () {
  "use strict";

  var TEAL = "#2b6a6a";
  var CLAY = "#c25a3c";
  var INK = "#17161a";
  var MUTE = "#8b877f";
  var RULE = "#e5e1d7";
  var LIGHT = "#c9c4b6";

  var reduce = window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function setup(id, h) {
    var c = document.getElementById(id);
    if (!c) return null;
    var dpr = window.devicePixelRatio || 1;
    function resize() {
      var w = c.clientWidth || c.parentNode.clientWidth || 640;
      c.width = Math.round(w * dpr);
      c.height = Math.round(h * dpr);
      c.style.height = h + "px";
      var g = c.getContext("2d");
      g.setTransform(dpr, 0, 0, dpr, 0, 0);
      c._w = w; c._h = h;
      return g;
    }
    return { c: c, resize: resize, h: h };
  }

  function loop(fn, staticT) {
    // Under reduced motion draw one frame chosen to be the most informative
    // state of that particular animation, not an arbitrary mid-transition.
    if (reduce) { fn(staticT); return; }
    function tick(ms) { fn((ms % 9000) / 9000); requestAnimationFrame(tick); }
    requestAnimationFrame(tick);
  }

  /* ---------------------------------------------------------------- points */

  // Deterministic cloud: two gaussian blobs plus a thin bridge, in unit space.
  function makeCloud(n) {
    var pts = [];
    function rnd(seed) { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
    var s = 7;
    function gauss() {
      var u = Math.max(rnd(s = (s + 1) | 0), 1e-6), v = rnd(s = (s + 1) | 0);
      return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    }
    for (var i = 0; i < n; i++) {
      var left = i % 2 === 0;
      var cx = left ? 0.34 : 0.66, cy = left ? 0.52 : 0.46;
      pts.push({
        x: cx + gauss() * 0.085,
        y: cy + gauss() * 0.11,
        left: left,
        ph: rnd(s = (s + 1) | 0) * 6.283
      });
    }
    return pts;
  }

  var cloud = makeCloud(150);

  var P = setup("c-points", 330);
  if (P) {
    var drawPoints = function (t) {
      var g = P.resize();
      var w = P.c._w, h = P.c._h;
      g.clearRect(0, 0, w, h);
      g.fillStyle = "#fff"; g.fillRect(0, 0, w, h);

      // Phase weights, each a smooth ramp.
      var ramp = function (a, b) { return Math.max(0, Math.min(1, (t - a) / (b - a))); };
      var wPts = 1 - ramp(0.30, 0.46);            // points fade as density arrives
      var wDen = ramp(0.24, 0.50) * (1 - ramp(0.78, 0.94));
      var wClu = ramp(0.62, 0.80);                 // clusters colour in
      var intro = ramp(0.0, 0.16);

      var pad = 26, W = w - pad * 2, H = h - pad * 2;
      var X = function (u) { return pad + u * W; };
      var Y = function (v) { return pad + v * H; };

      // density layer: additive soft blobs, drawn faintly
      if (wDen > 0.01) {
        g.globalAlpha = 0.16 * wDen;
        for (var i = 0; i < cloud.length; i++) {
          var p = cloud[i], R = 44;
          var grd = g.createRadialGradient(X(p.x), Y(p.y), 0, X(p.x), Y(p.y), R);
          grd.addColorStop(0, TEAL);
          grd.addColorStop(1, "rgba(43,106,106,0)");
          g.fillStyle = grd;
          g.beginPath(); g.arc(X(p.x), Y(p.y), R, 0, 6.2832); g.fill();
        }
        g.globalAlpha = 1;
      }

      // points
      for (var j = 0; j < cloud.length; j++) {
        var q = cloud[j];
        var jitter = reduce ? 0 : Math.sin(t * 12 + q.ph) * 1.4 * intro;
        var col;
        if (wClu > 0.5) col = q.left ? TEAL : CLAY;
        else col = "#6f6b62";
        g.globalAlpha = (0.85 * wPts + 0.9 * wClu * 0.35) * intro;
        g.fillStyle = col;
        g.beginPath();
        g.arc(X(q.x) + jitter, Y(q.y), wClu > 0.5 ? 3.4 : 2.7, 0, 6.2832);
        g.fill();
      }
      g.globalAlpha = 1;

      // cluster labels once they exist
      if (wClu > 0.6) {
        g.globalAlpha = (wClu - 0.6) / 0.4;
        g.font = "600 12px " + mono();
        g.textAlign = "center";
        g.fillStyle = TEAL; g.fillText("cluster A", X(0.34), Y(0.52) + 92);
        g.fillStyle = CLAY; g.fillText("cluster B", X(0.66), Y(0.46) + 92);
        g.globalAlpha = 1;
      }

      // stage caption, bottom-left
      g.font = "500 12px " + mono();
      g.textAlign = "left";
      g.fillStyle = MUTE;
      var cap = t < 0.26 ? "1 · points" : (t < 0.66 ? "2 · density" : "3 · clusters");
      g.fillText(cap, pad, h - 12);
    };
    loop(drawPoints, 0.95);   // static frame: clusters coloured
    window.addEventListener("resize", function () { drawPoints(0.95); });
  }

  /* ----------------------------------------------------------------- merge */

  // A one-dimensional density with two peaks and a valley between them.
  function f(x) {
    var g = function (c, s) { return Math.exp(-0.5 * Math.pow((x - c) / s, 2)); };
    return g(0.34, 0.095) + g(0.66, 0.095) + 0.06;
  }

  var M = setup("c-merge", 330);
  if (M) {
    // Valley height: sample the midpoint.
    var VALLEY = f(0.5);

    var drawMerge = function (t) {
      var g = M.resize();
      var w = M.c._w, h = M.c._h;
      g.clearRect(0, 0, w, h);
      g.fillStyle = "#fff"; g.fillRect(0, 0, w, h);

      var pad = 34, top = 26, base = h - 40;
      var W = w - pad * 2;
      var X = function (u) { return pad + u * W; };
      var peak = 1.13;                                  // a little headroom
      var Y = function (v) { return base - (v / peak) * (base - top); };

      // level sweeps from above the peaks to below the valley, then holds.
      var sweep = Math.min(1, t / 0.72);
      var eased = 1 - Math.pow(1 - sweep, 2);
      // Start just below the peaks (two caps) and end below the valley (one
      // region). Starting above the peaks would show an empty frame first.
      var level = 1.00 + (0.30 - 1.00) * eased;

      // filled region where f >= level, coloured by connected component
      var N = 240, runs = [], cur = null;
      for (var i = 0; i <= N; i++) {
        var x = i / N, above = f(x) >= level;
        if (above) { if (!cur) { cur = { a: x, b: x }; runs.push(cur); } else cur.b = x; }
        else cur = null;
      }

      // area under the curve, faint
      g.beginPath();
      g.moveTo(X(0), base);
      for (var k = 0; k <= N; k++) g.lineTo(X(k / N), Y(f(k / N)));
      g.lineTo(X(1), base); g.closePath();
      g.fillStyle = "#f6f4ee"; g.fill();

      // the above-level region, coloured per component
      for (var r = 0; r < runs.length; r++) {
        var run = runs[r];
        var col = runs.length > 1 ? (r === 0 ? TEAL : CLAY) : CLAY;
        g.globalAlpha = 0.22;
        g.beginPath();
        g.moveTo(X(run.a), base);
        for (var m = Math.ceil(run.a * N); m <= Math.floor(run.b * N); m++)
          g.lineTo(X(m / N), Y(f(m / N)));
        g.lineTo(X(run.b), base); g.closePath();
        g.fillStyle = col; g.fill();
        g.globalAlpha = 1;
      }

      // the curve itself
      g.beginPath();
      for (var k2 = 0; k2 <= N; k2++) {
        var px = X(k2 / N), py = Y(f(k2 / N));
        if (k2 === 0) g.moveTo(px, py); else g.lineTo(px, py);
      }
      g.strokeStyle = INK; g.lineWidth = 2; g.stroke();

      // the sweeping level
      g.beginPath(); g.moveTo(pad - 8, Y(level)); g.lineTo(w - pad + 8, Y(level));
      g.setLineDash([5, 5]); g.strokeStyle = "#9a958a"; g.lineWidth = 1.5; g.stroke();
      g.setLineDash([]);

      // valley marker
      g.beginPath(); g.moveTo(X(0.5), Y(VALLEY)); g.lineTo(X(0.5), base);
      g.strokeStyle = RULE; g.lineWidth = 1; g.stroke();

      // state caption
      var two = runs.length > 1;
      g.font = "600 12px " + mono();
      g.textAlign = "left";
      g.fillStyle = two ? TEAL : CLAY;
      g.fillText(runs.length === 0 ? "0 clusters" : (two ? "2 clusters" : "1 cluster"),
                 pad, h - 14);
      g.textAlign = "right";
      g.fillStyle = MUTE;
      g.fillText("level " + level.toFixed(2), w - pad, h - 14);

      g.textAlign = "center";
      g.fillStyle = MUTE;
      g.font = "500 11px " + mono();
      g.fillText("valley depth = persistence", X(0.5), h - 14);
    };
    loop(drawMerge, 0.30);    // static frame: level just above the valley
    window.addEventListener("resize", function () { drawMerge(0.30); });
  }

  function mono() {
    return 'ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace';
  }
})();
