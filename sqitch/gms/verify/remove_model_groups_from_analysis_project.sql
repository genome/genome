-- Verify remove_model_groups_from_analysis_project

DO $$
BEGIN

IF EXISTS(SELECT * from information_schema.columns
    WHERE table_schema = 'config'
    AND table_name = 'analysis_project'
    AND column_name = 'model_group_id') THEN
    RAISE EXCEPTION 'model_group_id still exists!';
END IF;

END $$
