-- Revert subject.subject_attribute.index_attribute_label_subject_id

BEGIN;

DROP INDEX subject.idx_s_sa_al_si;

COMMIT;
