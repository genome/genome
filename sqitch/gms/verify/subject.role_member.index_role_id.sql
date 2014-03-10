-- Verify subject.role_member.index_role_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_s_rm_ri';

ROLLBACK;
