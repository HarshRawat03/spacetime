# Spacetime

Interactive calculators for every chapter of O. Levenspiel, *Chemical Reaction Engineering*, 3rd ed. (Wiley, 1999): 67 tools covering chapters 1–30 and the appendix.

Each tool takes the chapter's ideas as far as they go numerically. Rate data analysis by every method in Chapter 3, ideal reactors with any kinetics, multiple reactions, non-isothermal design and staging, RTD and the flow models, pore diffusion with heat effects and egg-shell pellets, deactivation, gas–liquid reactions with the exact film-theory enhancement factor, shrinking-core design for size mixtures, and enzyme and microbial fermentation. The worked examples and many end-of-chapter problems load as presets, and the tool compares its answer with the printed one.

A project by Harsh Rawat, Chemical Engineering, IIT Jammu. Not affiliated with the author or the publisher.

## Build

```
npm install
node build.js        # writes dist/index.html, one self-contained file
```

Everything (fonts, KaTeX styles, pre-rendered equations, all scripts) is inlined, so `dist/index.html` can be opened from disk or served from any static host. Nothing is sent anywhere: all computation runs in the browser, and inputs are saved in the browser's local storage only.

## Tests

```
node test/t_num.js      # the numerical library: special functions, quadrature, ODEs, root finding
node test/run_tools.js  # every preset of every tool; pass a tag prefix to filter, e.g. "C-18"
```

`run_tools.js` reports three kinds of checks for each preset: the book's printed answer (`book`), the book's own formula evaluated exactly (`expect`), and values worked out independently where the book gives none (`check`).

## Layout

```
src/core/   numerical library (num.js), kinetics and networks, thermodynamics, RTD, pellet
            diffusion, plotting, P&ID drawing, the light math markup, the UI shell
src/tools/  one file per part or chapter: p1-* (Part I) … p5-* (Part V), pA-31 (appendix)
src/css/    the stylesheet
build.js    bundles everything into dist/index.html and checks every equation with KaTeX
```

A tool is a plain object passed to `ST.tool({...})`: inputs, presets, equations, a method note, and a `run(v)` function that returns key results, warnings, plots, tables and an optional process sketch.
