-- Deploy result.user.user_id
-- requires: result_user

BEGIN;

CREATE INDEX sru_uid_i on result."user" using btree (user_id);

COMMIT;
