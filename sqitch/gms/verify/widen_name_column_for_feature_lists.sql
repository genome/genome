-- Verify widen_name_column_for_feature_lists

BEGIN;

  SELECT name FROM model.feature_list WHERE FALSE;

ROLLBACK;
