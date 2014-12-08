-- Verify make_model_group_id_nullable
DO $$
BEGIN

  IF EXISTS(SELECT * from information_schema.columns
      WHERE table_schema = 'config'
      AND table_name = 'analysis_project'
      AND column_name = 'model_group_id'
      AND is_nullable = 'NO') THEN
      RAISE EXCEPTION 'model_group_id is not nullable';
  END IF;

END $$;
