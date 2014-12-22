-- Verify remove_params_and_inputs_from_software_results

DO $$
BEGIN

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'result'
    AND table_name = 'software_result'
    AND column_name = 'params_id') THEN
    RAISE EXCEPTION 'params_id still exists!';
END IF;

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'result'
    AND table_name = 'software_result'
    AND column_name = 'inputs_id') THEN
    RAISE EXCEPTION 'inputs still exists!';
END IF;

END $$;
