-- Revert subject.role_member.index_user_email

BEGIN;

DROP INDEX subject.idx_s_rm_ue;

COMMIT;
