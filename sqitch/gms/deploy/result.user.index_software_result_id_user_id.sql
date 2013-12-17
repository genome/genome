-- Deploy result.user.software_result_id_user_id
-- requires: result_user

BEGIN;

CREATE INDEX user_result_id_user_id_index on result."user" using btree (software_result_id, user_id);

COMMIT;
