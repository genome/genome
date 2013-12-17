-- Deploy model.model.genome_model_id_subject_class_name
-- requires: model_model

BEGIN;

CREATE INDEX m_m_id_subject_class_name_index on model.model using btree (genome_model_id, subject_class_name);

COMMIT;
