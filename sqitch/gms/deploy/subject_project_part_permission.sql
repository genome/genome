-- Deploy subject_project_part_permission
-- requires: subject_project_part

BEGIN;

REVOKE ALL ON TABLE subject.project_part FROM PUBLIC;
REVOKE ALL ON TABLE subject.project_part FROM genome;
GRANT ALL ON TABLE subject.project_part TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.project_part TO "gms-user";

COMMIT;
