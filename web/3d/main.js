/* Entry point of web/vacancy3d.js: the two 3D figures of the
   vacancy-distribution workspace. The page hands each figure the numbers
   it read from the committed results; nothing is computed here. */
import { webglAvailable } from './common.js';
import { createParticle, DOT_UMOL } from './particle.js';
import { createSlab, NX, NY, LAYERS } from './slab.js';

export const available = webglAvailable;
export const particle = createParticle;
export const slab = createSlab;
export const PATCH = { NX, NY, LAYERS, DOT_UMOL };
