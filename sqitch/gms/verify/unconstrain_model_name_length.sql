-- Verify unconstrain_model_name_length

BEGIN;

  SELECT 1/(character_maximum_length IS NULL)::int b
  FROM information_schema.columns
  WHERE table_schema = 'model' AND table_name = 'model' AND column_name = 'name';

ROLLBACK;
