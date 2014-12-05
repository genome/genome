-- Verify remove_width_constraint_for_result_users


BEGIN;

  SELECT 1/count(*) from information_schema.columns
      WHERE table_schema = 'result'
      AND table_name = 'user'
      AND column_name = 'label'
      AND data_type = 'text';

  SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'result_user_label_idx';

ROLLBACK;

