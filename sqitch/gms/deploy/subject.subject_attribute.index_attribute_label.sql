-- Deploy subject.subject_attribute.attribute_label
-- requires: subject_subject_attribute

BEGIN;

CREATE INDEX subject_attribute_label_index on subject.subject_attribute using btree (attribute_label);

COMMIT;
