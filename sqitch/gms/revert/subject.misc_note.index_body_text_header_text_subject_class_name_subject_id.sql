-- Revert subject.misc_note.index_body_text_header_text_subject_class_name_subject_id

BEGIN;

DROP INDEX subject.idx_s_mn_bt_ht_scn_si;

COMMIT;
