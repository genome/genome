-- Revert model.model.index_subject_class_name_subject_id

BEGIN;

DROP INDEX model.model_subject_index;

COMMIT;
