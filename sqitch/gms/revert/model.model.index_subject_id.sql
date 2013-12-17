-- Revert model.model.index_subject_id

BEGIN;

DROP INDEX model.model_subject_id_index;

COMMIT;
