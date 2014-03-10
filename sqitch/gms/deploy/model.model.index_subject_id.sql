-- Deploy model.model.subject_id
-- requires: model_model

BEGIN;

CREATE INDEX model_subject_id_index on model.model using btree (subject_id);

COMMIT;
