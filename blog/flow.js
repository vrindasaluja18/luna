/* LUNA — data-flow animation.
   Dependency-free canvas. Data particles travel left to right through the
   pipeline. Their speed is scaled by the *resident bytes* of the stage they are
   crossing, so they bunch up at the Morse topology stage — the 43 MB transient
   is the bottleneck, and the animation is meant to make that legible rather
   than decorative. Clicking a stage row below emphasises its node. */

(function () {
  "use strict";

  var canvas = document.getElementById("flowCanvas");
  if (!canvas) return;
  var ctx = canvas.getContext("2d");

  var reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  // weight drives node size and particle speed. topology = 43 MB is the peak.
  var STAGES = [
    { k: "sim",       nm: "sim",       bytes: "O(1)",       weight: 4   },
    { k: "frontend",  nm: "frontend",  bytes: "O(feat)",    weight: 8   },
    { k: "backend",   nm: "backend",   bytes: "O(kf)",      weight: 14  },
    { k: "mapping",   nm: "mapping",   bytes: "17.3 MB",    weight: 26  },
    { k: "topology",  nm: "topology",  bytes: "~43 MB",     weight: 100 },
    { k: "planning",  nm: "planning",  bytes: "O(nodes)",   weight: 10  },
    { k: "control",   nm: "control",   bytes: "O(1)",       weight: 6   }
  ];

  var C = {
    ink:   "#0a0c0f",
    rule:  "#232a33",
    dim:   "#aab4c0",
    paper: "#e9edf2",
    amber: "#ffb454",
    trace: "#66e0ff",
    hot:   "#ffb454"
  };

  var dpr = Math.min(window.devicePixelRatio || 1, 2);
  var W = 1100, H = 300;
  var nodes = [];
  var particles = [];
  var active = null;      // stage key currently emphasised
  var hover = null;

  var NPART = 64;

  function layout() {
    var padL = 78, padR = 78, padT = 118;
    var n = STAGES.length;
    var span = W - padL - padR;
    nodes = STAGES.map(function (s, i) {
      var x = padL + (span * i) / (n - 1);
      return {
        k: s.k, nm: s.nm, bytes: s.bytes, weight: s.weight,
        x: x, y: padT, r: 5 + Math.sqrt(s.weight) * 1.05
      };
    });
    // cumulative path lengths for particle positioning
    var acc = 0;
    nodes[0].s0 = 0;
    for (var i = 1; i < nodes.length; i++) {
      acc += Math.hypot(nodes[i].x - nodes[i - 1].x, nodes[i].y - nodes[i - 1].y);
      nodes[i].s0 = acc;
    }
    nodes.total = acc;
  }

  function pointAt(s) {
    var total = nodes.total;
    s = Math.max(0, Math.min(total, s));
    for (var i = 1; i < nodes.length; i++) {
      if (s <= nodes[i].s0) {
        var a = nodes[i - 1], b = nodes[i];
        var seg = nodes[i].s0 - a.s0;
        var t = seg > 0 ? (s - a.s0) / seg : 0;
        return { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t, i: i };
      }
    }
    var last = nodes[nodes.length - 1];
    return { x: last.x, y: last.y, i: nodes.length - 1 };
  }

  // speed multiplier: heavy stages slow the flow. normalised so 43 MB -> 0.22
  function speedAt(s) {
    var p = pointAt(s);
    var near = 0;
    for (var i = 0; i < nodes.length; i++) {
      var d = Math.abs(p.x - nodes[i].x);
      var reach = nodes[i].r + 26;
      if (d < reach) near = Math.max(near, nodes[i].weight / 100 * (1 - d / reach));
    }
    return 1 - 0.78 * near;
  }

  function resetParticle(p) {
    p.s = Math.random() * nodes.total;
    p.base = 42 + Math.random() * 46;   // px per second
    p.alpha = 0.35 + Math.random() * 0.5;
  }

  function initParticles() {
    particles = [];
    for (var i = 0; i < NPART; i++) {
      var p = {};
      resetParticle(p);
      p.s = (i / NPART) * nodes.total;   // even spread on first frame
      particles.push(p);
    }
  }

  function drawNode(n) {
    var emphasised = active === n.k || hover === n.k;
    var isHot = n.k === "topology";
    var stroke = emphasised ? C.amber : (isHot ? C.trace : "#3a4657");
    var fill = emphasised ? "rgba(255,180,84,.16)"
             : (isHot ? "rgba(102,224,255,.12)" : "rgba(255,255,255,.03)");

    // square node (blueprint, not a bubble)
    var r = n.r;
    ctx.fillStyle = fill;
    ctx.strokeStyle = stroke;
    ctx.lineWidth = emphasised ? 2 : 1.2;
    ctx.beginPath();
    ctx.rect(n.x - r, n.y - r, r * 2, r * 2);
    ctx.fill();
    ctx.stroke();

    // resident-bytes bar under the node
    var bw = Math.max(4, n.r * 2);
    var bh = 3;
    var frac = n.weight / 100;
    ctx.fillStyle = "#232a33";
    ctx.fillRect(n.x - bw / 2, n.y + r + 8, bw, bh);
    ctx.fillStyle = isHot ? C.trace : C.amber;
    ctx.fillRect(n.x - bw / 2, n.y + r + 8, bw * frac, bh);

    // labels
    ctx.textAlign = "center";
    ctx.font = "11px ui-monospace, Menlo, monospace";
    ctx.fillStyle = emphasised ? C.amber : C.paper;
    ctx.fillText(n.nm, n.x, n.y + r + 26);
    ctx.font = "10px ui-monospace, Menlo, monospace";
    ctx.fillStyle = isHot ? C.trace : C.dim;
    ctx.fillText(n.bytes, n.x, n.y + r + 41);
  }

  function draw(t) {
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = C.ink;
    ctx.fillRect(0, 0, W, H);

    // faint vertical guides
    ctx.strokeStyle = C.rule;
    ctx.lineWidth = 1;
    ctx.beginPath();
    nodes.forEach(function (n) {
      ctx.moveTo(n.x, 34);
      ctx.lineTo(n.x, H - 42);
    });
    ctx.stroke();

    // the rail
    ctx.beginPath();
    ctx.moveTo(nodes[0].x, nodes[0].y);
    for (var i = 1; i < nodes.length; i++) ctx.lineTo(nodes[i].x, nodes[i].y);
    ctx.strokeStyle = "#2a323d";
    ctx.lineWidth = 1.5;
    ctx.stroke();

    // particles
    for (var j = 0; j < particles.length; j++) {
      var p = particles[j];
      var pos = pointAt(p.s);
      var tail = pointAt(p.s - 16);
      ctx.strokeStyle = "rgba(102,224,255," + p.alpha.toFixed(2) + ")";
      ctx.lineWidth = 1.6;
      ctx.beginPath();
      ctx.moveTo(tail.x, tail.y);
      ctx.lineTo(pos.x, pos.y);
      ctx.stroke();
    }

    nodes.forEach(drawNode);

    // caption
    ctx.textAlign = "left";
    ctx.font = "10px ui-monospace, Menlo, monospace";
    ctx.fillStyle = C.dim;
    ctx.fillText("data flow →", 14, 20);
    ctx.textAlign = "right";
    ctx.fillStyle = C.trace;
    ctx.fillText("particles slow where bytes are resident", W - 14, 20);

    // header
    ctx.textAlign = "left";
    ctx.font = "11px ui-monospace, Menlo, monospace";
    ctx.fillStyle = C.paper;
    ctx.fillText("resident bytes per stage", 14, H - 16);
    ctx.fillStyle = C.hot;
    ctx.fillText("● topology ≈ 43 MB — the largest transient in the pipeline", 14, H - 2);
  }

  var last = 0;
  function tick(ts) {
    if (!last) last = ts;
    var dt = Math.min(0.05, (ts - last) / 1000);
    last = ts;

    for (var i = 0; i < particles.length; i++) {
      var p = particles[i];
      p.s += p.base * speedAt(p.s) * dt;
      if (p.s >= nodes.total) { p.s = 0; resetParticle(p); p.s = 0; }
    }
    draw(ts);
    requestAnimationFrame(tick);
  }

  function resize() {
    var cssW = canvas.clientWidth || W;
    var scale = cssW / W;
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    canvas.style.height = (H * scale) + "px";
  }

  // interaction: stage rows emphasise their node
  var rows = document.querySelectorAll(".stage");
  rows.forEach(function (row) {
    var k = row.getAttribute("data-k");
    row.addEventListener("mouseenter", function () { hover = k; });
    row.addEventListener("mouseleave", function () { if (hover === k) hover = null; });
    row.addEventListener("click", function () {
      active = (active === k) ? null : k;
      rows.forEach(function (r) { r.classList.toggle("active", r.getAttribute("data-k") === active); });
    });
  });

  var yr = document.getElementById("yr");
  if (yr) yr.textContent = new Date().getFullYear();

  layout();
  initParticles();
  resize();
  window.addEventListener("resize", function () { resize(); draw(0); });

  if (reduced) {
    // one static frame, particles spread along the path
    draw(0);
  } else {
    requestAnimationFrame(tick);
  }
})();
