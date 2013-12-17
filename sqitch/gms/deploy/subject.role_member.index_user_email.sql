-- Deploy subject.role_member.user_email
-- requires: subject_role_member

BEGIN;

CREATE INDEX idx_s_rm_ue on subject.role_member using btree (user_email);

COMMIT;
