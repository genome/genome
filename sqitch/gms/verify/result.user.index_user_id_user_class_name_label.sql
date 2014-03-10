-- Verify result.user.index_user_id_user_class_name_label

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'user_id_name_label_index';

ROLLBACK;
