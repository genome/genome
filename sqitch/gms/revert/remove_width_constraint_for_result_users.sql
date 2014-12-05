-- Revert remove_width_constraint_for_result_users

BEGIN;

  ALTER TABLE result.user ALTER COLUMN label TYPE character varying(100);
  DROP INDEX result.result_user_label_idx;

COMMIT;
