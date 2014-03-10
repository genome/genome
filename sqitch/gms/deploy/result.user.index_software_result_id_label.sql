-- Deploy result.user.software_result_id_label
-- requires: result_user

BEGIN;

CREATE INDEX user_result_label_index on result."user" using btree (software_result_id, label);

COMMIT;
