-- Revert config.tag_subject_mapping

BEGIN;

DROP TABLE config.tag_subject_mapping;

COMMIT;
