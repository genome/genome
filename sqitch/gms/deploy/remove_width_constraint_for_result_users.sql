-- Deploy remove_width_constraint_for_result_users
-- requires: result_user

BEGIN;

  ALTER TABLE result.user ALTER COLUMN label TYPE text;

COMMIT;

CREATE INDEX CONCURRENTLY result_user_label_idx ON result.user (label);

