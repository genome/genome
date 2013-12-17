-- Deploy result.user.user_id_user_class_name_label
-- requires: result_user

BEGIN;

CREATE INDEX user_id_name_label_index on result."user" using btree (user_id, user_class_name, label);

COMMIT;
