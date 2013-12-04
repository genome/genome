-- Deploy result.user.software_result_id
-- requires: result_user

BEGIN;

CREATE INDEX user_software_result_id_idx on result."user" using btree (software_result_id);

COMMIT;
