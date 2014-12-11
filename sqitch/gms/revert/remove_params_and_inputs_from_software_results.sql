-- Revert remove_params_and_inputs_from_software_results

BEGIN;

  ALTER TABLE result.software_result ADD COLUMN inputs_id character varying(4000);
  ALTER TABLE result.software_result ADD COLUMN params_id character varying(4000);

COMMIT;
