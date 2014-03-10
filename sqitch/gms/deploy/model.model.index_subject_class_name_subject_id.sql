-- Deploy model.model.subject_class_name_subject_id
-- requires: model_model

BEGIN;

CREATE INDEX model_subject_index on model.model using btree (subject_class_name, subject_id);

COMMIT;
