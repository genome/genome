-- Revert subject.role_member.index_role_id

BEGIN;

DROP INDEX subject.idx_s_rm_ri;

COMMIT;
