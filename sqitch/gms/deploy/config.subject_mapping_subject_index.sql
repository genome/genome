-- Deploy config.subject_mapping_subject_index
-- requires: config.subject_mapping

BEGIN;

CREATE INDEX subject_mapping_subject_subject_mapping_id_idx ON config.subject_mapping_subject USING btree (subject_mapping_id);

COMMIT;
