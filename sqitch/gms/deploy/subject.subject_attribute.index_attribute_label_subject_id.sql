-- Deploy subject.subject_attribute.attribute_label_subject_id
-- requires: subject_subject_attribute

BEGIN;

CREATE INDEX idx_s_sa_al_si on subject.subject_attribute using btree (attribute_label, subject_id);

COMMIT;
