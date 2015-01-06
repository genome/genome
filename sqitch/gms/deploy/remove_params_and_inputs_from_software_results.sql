-- Deploy remove_params_and_inputs_from_software_results
-- requires: result_software_result

BEGIN;

  ALTER TABLE result.software_result DROP COLUMN inputs_id;
  ALTER TABLE result.software_result DROP COLUMN params_id;

COMMIT;
