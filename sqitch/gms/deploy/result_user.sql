-- Deploy result_user
-- requires: result_software_result

BEGIN;

CREATE TABLE IF NOT EXISTS result."user" (
    id character varying(32) NOT NULL,
    software_result_id character varying(32) NOT NULL,
    user_id character varying(255) NOT NULL,
    user_class_name character varying(255),
    label character varying(100),
    active boolean DEFAULT true NOT NULL,
    CONSTRAINT user_pkey PRIMARY KEY (id),
    CONSTRAINT fk_r_u_srid FOREIGN KEY (software_result_id) REFERENCES result.software_result(id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
