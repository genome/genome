-- Verify config.subject_mapping

BEGIN;

SELECT id,created_at,updated_at,created_by,analysis_project_id FROM config.subject_mapping WHERE FALSE;
SELECT id,created_at,updated_at,created_by,subject_mapping_id,key,value FROM config.subject_mapping_input WHERE FALSE;
SELECT id,created_at,updated_at,created_by,subject_mapping_id,label,subject_id FROM config.subject_mapping_subject WHERE FALSE;

ROLLBACK;
