-- Revert subject.subject_attribute.index_attribute_value_nomenclature

BEGIN;

DROP INDEX subject.idx_s_sa_av_n;

COMMIT;
