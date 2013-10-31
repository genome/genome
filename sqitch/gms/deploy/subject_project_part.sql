-- Deploy subject_project_part
-- requires: subject_project

BEGIN;

CREATE TABLE IF NOT EXISTS subject.project_part (
    id character varying(64) NOT NULL,
    project_id character varying(64) NOT NULL,
    part_class_name character varying(256) NOT NULL,
    part_id character varying(64) NOT NULL,
    label character varying(100),
    role character varying(100),
    CONSTRAINT project_part_pkey PRIMARY KEY (id),
    CONSTRAINT project_part_project_id_part_class_name_part_id_role_key
        UNIQUE (project_id, part_class_name, part_id, role),
    CONSTRAINT project_part_project_id_fkey FOREIGN KEY (project_id) REFERENCES subject.project(id)
);

COMMIT;
