# ExposureBounds

ExposureBounds certifies what partial measurements can establish about a delivered temporal signal. It handles bounded point readings, overlapping measurement windows and uncertain clock offsets. The verifier checks rational certificates of feasibility, exclusion and response sign over a continuous rate band.

## Reproduce

Use Python 3.14. From the repository root:

```sh
python -m pip install -r reproduction/requirements.txt
python -B reproduce_interface.py
cd reproduction
python -B -m unittest -q test_frontiers
python -B compare_band_bounds.py
```

The first script reproduces the exact per-arm interface example. The remaining commands verify the formal cases and compare global and segment-resolved spectral bounds. For all six cases, see [Guided reproduction](reproduction/GUIDED_REPRODUCTION.md). The archived input and output objects are formal examples.

## Castellano

ExposureBounds certifica qué se puede concluir sobre una señal temporal entregada a partir de mediciones parciales. Admite lecturas puntuales acotadas, ventanas solapadas e incertidumbre de reloj. El verificador comprueba certificados racionales de factibilidad, exclusión y signo de la respuesta sobre una banda continua de tasas.

Los comandos anteriores reproducen la interfaz exacta por brazo, verifican los casos formales y comparan las cotas espectrales globales y por segmento. La [guía bilingüe](reproduction/GUIDED_REPRODUCTION.md) explica las seis decisiones.

## Citation / Cita

Facundo Firmenich, *ExposureBounds: software and formal certificates for temporal exposure under acquisition and clock uncertainty*. Cite the archived version used in the analysis. / Citá la versión archivada que utilices en el análisis.

[ChronoWitness](https://github.com/FacundoFirmenich/ChronoWitness) uses per-arm delivery bounds in its prospective observation model. / ChronoWitness usa las cotas de entrega por brazo en su modelo prospectivo de observación.
