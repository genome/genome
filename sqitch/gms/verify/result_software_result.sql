-- Verify result_software_result

BEGIN;

SELECT id, class_name, version, inputs_id, params_id, outputs_path, lookup_hash
FROM result.software_result
WHERE FALSE;

ROLLBACK;
