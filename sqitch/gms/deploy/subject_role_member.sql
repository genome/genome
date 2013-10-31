-- Deploy subject_role_member
-- requires: subject_role

BEGIN;

CREATE TABLE IF NOT EXISTS subject.role_member (
    user_email character varying(255) NOT NULL,
    role_id character varying(32) NOT NULL,
    CONSTRAINT role_member_pkey PRIMARY KEY (user_email, role_id),
    CONSTRAINT role_member_role_id_fkey FOREIGN KEY (role_id) REFERENCES subject.role(id),
    CONSTRAINT role_member_user_email_fkey FOREIGN KEY (user_email) REFERENCES subject."user"(email)
);

COMMIT;
