-- Verify subject.role_member.index_user_email

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_s_rm_ue';

ROLLBACK;
