-- Verify config.tag_subject_mapping

BEGIN;

SELECT *
FROM config.tag_subject_mapping
WHERE FALSE;

ROLLBACK;
