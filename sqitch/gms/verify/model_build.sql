-- Verify model_build

BEGIN;

SELECT build_id, data_directory, model_id, software_revision, subclass_name
FROM model.build
WHERE FALSE;

ROLLBACK;
