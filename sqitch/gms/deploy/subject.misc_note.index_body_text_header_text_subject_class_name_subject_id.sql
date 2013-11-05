-- Deploy subject.misc_note.body_text_header_text_subject_class_name_subject_id
-- requires: subject_misc_note

BEGIN;

CREATE INDEX idx_s_mn_bt_ht_scn_si on subject.misc_note using btree (body_text, header_text, subject_class_name, subject_id);

COMMIT;
