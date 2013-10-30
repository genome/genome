-- Verify model_feature_list

BEGIN;

SELECT id, name, source, format, reference_id, subject_id,
    file_id, file_content_hash,content_type, description
FROM model.feature_list
WHERE FALSE;

ROLLBACK;
