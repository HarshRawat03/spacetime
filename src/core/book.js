/* Structure of Levenspiel, Chemical Reaction Engineering, 3rd ed. (Wiley, 1999). */
(function (G) {
  'use strict';
  const ST = G.ST || (G.ST = {});
  ST.BOOK = {
    parts: [
      { n: '0', code: 'O', t: 'Overview', ch: [1] },
      { n: 'I', code: 'H', t: 'Homogeneous Reactions in Ideal Reactors', ch: [2, 3, 4, 5, 6, 7, 8, 9, 10] },
      { n: 'II', code: 'F', t: 'Flow Patterns, Contacting, and Non-Ideal Flow', ch: [11, 12, 13, 14, 15, 16] },
      { n: 'III', code: 'C', t: 'Reactions Catalyzed by Solids', ch: [17, 18, 19, 20, 21, 22] },
      { n: 'IV', code: 'N', t: 'Non-Catalytic Systems', ch: [23, 24, 25, 26] },
      { n: 'V', code: 'B', t: 'Biochemical Reaction Systems', ch: [27, 28, 29, 30] },
      { n: 'A', code: 'A', t: 'Appendix', ch: [31] }
    ],
    ch: {
      1: { t: 'Overview of Chemical Reaction Engineering', p: 1 },
      2: { t: 'Kinetics of Homogeneous Reactions', p: 13 },
      3: { t: 'Interpretation of Batch Reactor Data', p: 38 },
      4: { t: 'Introduction to Reactor Design', p: 83 },
      5: { t: 'Ideal Reactors for a Single Reaction', p: 90 },
      6: { t: 'Design for Single Reactions', p: 120 },
      7: { t: 'Design for Parallel Reactions', p: 152 },
      8: { t: 'Potpourri of Multiple Reactions', p: 170 },
      9: { t: 'Temperature and Pressure Effects', p: 207 },
      10: { t: 'Choosing the Right Kind of Reactor', p: 240 },
      11: { t: 'Basics of Non-Ideal Flow', p: 257 },
      12: { t: 'Compartment Models', p: 283 },
      13: { t: 'The Dispersion Model', p: 293 },
      14: { t: 'The Tanks-in-Series Model', p: 321 },
      15: { t: 'The Convection Model for Laminar Flow', p: 339 },
      16: { t: 'Earliness of Mixing, Segregation, and RTD', p: 350 },
      17: { t: 'Heterogeneous Reactions: Introduction', p: 369 },
      18: { t: 'Solid Catalyzed Reactions', p: 376 },
      19: { t: 'The Packed Bed Catalytic Reactor', p: 427 },
      20: { t: 'Reactors with Suspended Solid Catalyst, Fluidized Reactors of Various Types', p: 447 },
      21: { t: 'Deactivating Catalysts', p: 473 },
      22: { t: 'G/L Reactions on Solid Catalysts: Trickle Beds, Slurry Reactors, Three-Phase Fluidized Beds', p: 500 },
      23: { t: 'Fluid–Fluid Reactions: Kinetics', p: 523 },
      24: { t: 'Fluid–Fluid Reactors: Design', p: 540 },
      25: { t: 'Fluid–Particle Reactions: Kinetics', p: 566 },
      26: { t: 'Fluid–Particle Reactors: Design', p: 589 },
      27: { t: 'Enzyme Fermentation', p: 611 },
      28: { t: 'Microbial Fermentation: Introduction and Overall Picture', p: 623 },
      29: { t: 'Substrate-Limiting Microbial Fermentation', p: 630 },
      30: { t: 'Product-Limiting Microbial Fermentation', p: 645 },
      31: { t: 'Miscellany: Units, Magnitudes and Dimensionless Groups', p: 655 }
    }
  };
  ST.partOf = ch => ST.BOOK.parts.find(p => p.ch.indexOf(ch) >= 0);
  if (typeof module !== 'undefined') module.exports = ST;
})(typeof window !== 'undefined' ? window : globalThis);
