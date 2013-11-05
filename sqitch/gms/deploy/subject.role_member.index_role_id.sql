-- Deploy subject.role_member.role_id
-- requires: subject_role_member

BEGIN;

CREATE INDEX idx_s_rm_ri on subject.role_member using btree (role_id);

COMMIT;
