-- Deploy subject.user.username
-- requires: subject_user

BEGIN;

CREATE INDEX subject_user_username_index on subject."user" using btree (username);

COMMIT;
