# Guided reproduction | Reproducción guiada

## English

Start with the formal per-arm interface: run python -B reproduce_interface.py from the author package. Check that the two maximum input deviations are 1/8 and 1/16 and that their reporter-weighted sum is 3/(16e). The envelope figure visualises these intervals for the special case K=0.

From reproduction/, run python -B -m unittest -q test_frontiers and python -B compare_band_bounds.py. The latter compares global and segment-resolved certificates on the same spectral case. Compare the one-cell results before changing the grid.

Run python -B reproduce.py in a fresh unpacked copy and inspect results/FORMAL_RESULTS.json. Read its six cases in order: common-clock incompatibility, feasible boundary offsets, full exclusion of individual offsets, positive spectral band with a negative terminal tail, abstention under a small clock budget, and the conditional spectral mixture. For each, check input_parameters and the decision object together. The first three clock decisions are impossibility, feasibility and impossibility. The fifth case remains unresolved. The mixture additionally requires a declared coverage bound.

## Castellano

Empezá por la interfaz formal por brazo. Ejecutá python -B reproduce_interface.py desde el paquete de autor. Las desviaciones máximas son 1/8 y 1/16 y su suma ponderada por el reportero es 3/(16e). La figura de envolventes visualiza esos intervalos para K=0.

Desde reproduction/, ejecutá python -B -m unittest -q test_frontiers y python -B compare_band_bounds.py. El segundo comando compara los certificados globales y por segmento para el mismo caso espectral. Contrastá los resultados con una celda antes de cambiar la malla.

Ejecutá python -B reproduce.py en una copia recién descomprimida y abrí results/FORMAL_RESULTS.json. Leé los seis casos en orden: incompatibilidad con reloj común, offsets de frontera factibles, exclusión de toda la caja de offsets individuales, banda positiva con cola terminal negativa, abstención por presupuesto pequeño y mezcla espectral condicional. Contrastá input_parameters con la decisión en cada caso. Las tres primeras decisiones de reloj son imposibilidad, factibilidad e imposibilidad. El quinto caso queda sin resolver. La mezcla exige además una cota de cobertura declarada.
