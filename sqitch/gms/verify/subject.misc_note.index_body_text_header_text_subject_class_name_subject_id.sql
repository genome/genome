-- Verify subject.misc_note.index_body_text_header_text_subject_class_name_subject_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_s_mn_bt_ht_scn_si';

ROLLBACK;
