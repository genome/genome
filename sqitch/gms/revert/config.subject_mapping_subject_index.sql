-- Revert config.subject_mapping_subject_index

BEGIN;

DROP INDEX config.subject_mapping_subject_subject_mapping_id_idx;

COMMIT;
