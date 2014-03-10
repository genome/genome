-- Deploy subject.subject_attribute.attribute_value_nomenclature
-- requires: subject_subject_attribute

BEGIN;

CREATE INDEX idx_s_sa_av_n on subject.subject_attribute using btree (attribute_value, nomenclature);

COMMIT;
