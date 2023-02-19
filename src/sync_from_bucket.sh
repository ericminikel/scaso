#!/bin/bash
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_code/scaso_terra.R ./src/scaso_terra.R
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/figure_1.png ./display_items/figure-1.png
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/ana_cyno_cerebellum_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/ana_cyno_cortex_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/ana_mouse_cerebellum_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/ana_mouse_cortex_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/ana_mouse_thalamus_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mar_cyno_cerebellum_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mar_cyno_cortex_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mar_mouse_cerebellum_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mar_mouse_cortex_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mar_mouse_thalamus_full.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mouse_exc_clust.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/mouse_inh_clust.tsv ./data/analytical/
gsutil -m cp gs://fc-77539204-5aa7-4c5a-9286-e00ec3aa11b3/r_output/zeroes.tsv ./data/analytical/
