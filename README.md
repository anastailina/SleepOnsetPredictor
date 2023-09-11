
# MSc Project: Prediction of Sleep Onset using Feature-based Representation of EEG

This repository contains the code used for the completion of the MSc project by Anastasia Ilina. The project focuses on the prediction of sleep onset using feature-based representation of EEG through Random Forest and LSTM models on the MESA dataset.

## Table of Contents
1. [Technologies Used](#technologies-used)
2. [Matlab Component](#matlab-component)
3. [Python Component](#python-component)
4. [Folders Structure](#folders-structure)

---

## Technologies Used

- Python
- Matlab

---

## Matlab Component

The Matlab part of the project is responsible for:

- EEG and sleep scores extraction from MESA EDF files
- Pre-processing of EEG 
  - Butterworth filtering
  - Cardiac artefact correction
  - General artefact removal
- EEG epoching
- Feature extraction from EEG epochs

### Folders in Matlab Component

- `borrowed_work`: Contains work done prior to the start of this MSc project.
- `sleep_onset_analysis`: Contains additions carried out by the author during the course of this analysis.

---

## Python Component

The Python part of the project is comprised of Jupyter notebooks containing code for:

- Data exploration
- Data preprocessing for machine learning and deep learning tasks
- Model creation
- Training
- Hyperparameter tuning
- Final evaluation

## Folders in Python Component

- `final_models`: Contains final models selected for final tests on the test sets: non-stateful LSTM for classificationl, non-stateful LSTM for regression, Balanced Random Forest for Classification, along with the transfer learning exploration for regression task from final Balanced random forest and LSTM classifier models. In the Balanced Random Forest notebooks, the feature importance analysis for each timepoint prior to time to sleep onset is also implemented, but was not included in the final report. 
- `data_and_feature_analysis`: Contains the notebook for the exploration of the final feature-based EEG data.
- `model_experimentation`: Contains various experiments not as extensive as final models.
  - `on_final_data`: Contains experiments conducted on the data obtained from the final EEG pre-processing and feature extraction Matlab pipeline.
    - `covariate_shift`: Contains notebook summarizing the quantification of covariate shift between participants.
  - `on_unfixed_data`: Contains analysis on feature-based EEG representations with excessive general artefact rejection.
  - `early_experiments`: Contains early experiments with erroneous feature calculations.


---

